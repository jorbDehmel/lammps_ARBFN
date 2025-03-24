/**
 * @file controller.hpp
 * @brief Controller resources
 */

#pragma once

#include <boost/json/object.hpp>
#include <boost/json/src.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <mpi.h>
#include <thread>

/**
 * @brief A controller wherein every atom's fix is independent
 * of every other atom's. This is much more efficient than a
 * dependent controller if dependent control is not needed.
 * @param _single_atom_lambda The fix to call on every atom,
 * with the atom's data being the JSON first arg and the
 * resultant force deltas being saved in the second-fourth args.
 */
inline void independent_controller(
    std::function<void(const boost::json::object &, double &, double &, double &)>
        _single_atom_lambda,
    const uint64_t &_max_ms = 10000)
{
  MPI_Comm comm, junk_comm;
  MPI_Init(NULL, NULL);
  MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &junk_comm);
  MPI_Comm_split(MPI_COMM_WORLD, 56789, 0, &comm);

  std::cerr << __FILE__ << ":" << __LINE__ << "> "
            << "Started controller.\n"
            << std::flush;

  // For as long as there are connections left
  uint request_instance_counter = 0, requests = 0;
  uint num_registered = 0;
  uint64_t ms_since_update = 0;
  do {
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &status);

    if (flag) {
      ms_since_update = 0;
      char *const buffer = new char[status._ucount + 1];
      MPI_Recv(buffer, status._ucount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &status);
      buffer[status._ucount] = '\0';
      boost::json::object json = boost::json::parse(buffer).as_object();
      delete[] buffer;

      // Bookkeeping
      if (json["type"] == "register") {
        ++num_registered;
        const std::string raw = "{\"type\": \"ack\"}";
        MPI_Send(raw.c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
      } else if (json["type"] == "deregister") {
        assert(num_registered > 0);
        --num_registered;
      }

      // Data processing
      else if (json["type"] == "request") {
        // Periodically update user
        ++request_instance_counter;
        if (request_instance_counter % num_registered == 0) {
          request_instance_counter = 0;
          ++requests;

          if (requests % 100 == 0) { std::cerr << "Request #" << requests << "\n" << std::flush; }
        }

        // Determine fix to send back
        boost::json::array list;
        for (const auto &item : json["atoms"].as_array()) {
          double dfx = 0.0, dfy = 0.0, dfz = 0.0;

          // Processing here
          _single_atom_lambda(item.as_object(), dfx, dfy, dfz);

          boost::json::object fix;
          fix["dfx"] = dfx;
          fix["dfy"] = dfy;
          fix["dfz"] = dfz;
          list.push_back(fix);
        }

        // Properly format the response
        boost::json::object json_to_send;
        json_to_send["type"] = "response";
        json_to_send["atoms"] = list;
        std::stringstream s;
        s << json_to_send;
        const std::string raw = s.str();
        MPI_Send(s.str().c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
      }
    } else {
      // Delay
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      ms_since_update += 10;

      if (ms_since_update > _max_ms) { MPI_Abort(comm, 10); }
    }
  } while (num_registered != 0);

  std::cerr << __FILE__ << ":" << __LINE__ << "> "
            << "Halting controller\n"
            << std::flush;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&comm);
  MPI_Comm_free(&junk_comm);
  MPI_Finalize();
}

/**
 * @brief A prebuilt controller wherein no atom's force deltas
 * can be determined before all data has been reported. This is
 * less efficient than the independent route, but is necessary
 * in some cases (EG AI control).
 * @param _on_recv_all Called after all atoms have been
 * received. This provides all the data, the other lambda must
 * index into it. Called repeatedly (after some delay) until it
 * returns true. This makes it a good place to CHECK ON THE
 * PROGRESS OF expensive calls, but NOT a good place to wait on
 * them!
 * @param _single_atom Called for each atom after the other
 * lambda. Provides only an index, so you best hang onto the
 * bulk atom data provided in the previous callback.
 */
inline void dependent_controller(
    std::function<bool(const boost::json::array &)> _on_recv_all,
    std::function<void(const uint64_t &, double &, double &, double &)> _single_atom,
    const uint64_t &_max_ms = 10000)
{
  MPI_Comm comm, junk_comm;
  MPI_Init(NULL, NULL);
  MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &junk_comm);
  MPI_Comm_split(MPI_COMM_WORLD, 56789, 0, &comm);

  std::cerr << __FILE__ << ":" << __LINE__ << "> "
            << "Started controller.\n"
            << std::flush;

  // For as long as there are connections left
  uint requests = 0;
  uint num_registered = 0;

  // Bulk controlling: Maps worker rank to atom data
  // Cleared after every successful step
  std::map<int, boost::json::array> bulk_received;
  uint64_t ms_since_update = 0;

  do {
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &status);

    if (flag) {
      char *const buffer = new char[status._ucount + 1];
      MPI_Recv(buffer, status._ucount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &status);
      buffer[status._ucount] = '\0';
      boost::json::object json = boost::json::parse(buffer).as_object();
      delete[] buffer;

      // Bookkeeping
      if (json["type"] == "register") {
        ++num_registered;
        const std::string raw = "{\"type\": \"ack\"}";
        MPI_Send(raw.c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
      } else if (json["type"] == "deregister") {
        assert(num_registered > 0);
        --num_registered;
      }

      // Data processing
      else if (json["type"] == "request") {
        // Synchronization stuff
        bulk_received[status.MPI_SOURCE] = json.at("atoms").as_array();
        if (bulk_received.size() != num_registered) {
          // Send waiting packet and continue
          const std::string msg = "{\"type\": \"waiting\"}";
          MPI_Send(msg.c_str(), msg.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
          continue;
        }

        // Periodically update user
        ++requests;
        if ((requests / num_registered) % 1000 == 0) {
          std::cerr << "Request #" << (requests / num_registered) << "\n" << std::flush;
        }

        // Prepare list of all atoms
        boost::json::array list_to_send;
        for (const auto &p : bulk_received) {
          for (const auto &item : p.second) {
            // `item` is a single atom
            list_to_send.push_back(item.as_object());
          }
        }

        // Call first lambda until it returns true
        while (!_on_recv_all(list_to_send)) {
          for (const auto &p : bulk_received) {
            // Send waiting packet and continue
            const std::string msg = "{\"type\": \"waiting\"}";
            MPI_Send(msg.c_str(), msg.size(), MPI_CHAR, p.first, 0, comm);
          }
        }

        // Get atom info from second lambda and send to workers
        uint index = 0;
        for (const auto &p : bulk_received) {
          boost::json::array list;
          for (const auto &item : p.second) {
            double dfx = 0.0, dfy = 0.0, dfz = 0.0;

            // Processing here
            _single_atom(index, dfx, dfy, dfz);

            boost::json::object fix;
            fix["dfx"] = dfx;
            fix["dfy"] = dfy;
            fix["dfz"] = dfz;
            list.push_back(fix);

            ++index;
          }

          // Properly format the response
          boost::json::object json_to_send;
          json_to_send["type"] = "response";
          json_to_send["atoms"] = list;
          std::stringstream s;
          s << json_to_send;
          const std::string raw = s.str();
          MPI_Send(s.str().c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
        }

        bulk_received.clear();
      }
    } else {
      // Delay
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      ms_since_update += 10;

      if (ms_since_update > _max_ms) { MPI_Abort(comm, 10); }
    }
  } while (num_registered != 0);

  std::cerr << __FILE__ << ":" << __LINE__ << "> "
            << "Halting controller\n"
            << std::flush;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&comm);
  MPI_Comm_free(&junk_comm);
  MPI_Finalize();
}

/**
 * @brief ffield controller for more efficient special cases
 * @param _get_forces Maps array of atoms, flag indicating first
 * node, and position to forces (last argument).
 */
inline void ffield_controller(
    std::function<void(const boost::json::value &, const bool &, const double[3], double[3])>
        _get_forces)
{
  MPI_Comm comm, junk_comm;
  uintmax_t num_registered = 0;
  MPI_Init(NULL, NULL);
  MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &junk_comm);
  MPI_Comm_split(MPI_COMM_WORLD, 56789, 0, &comm);
  do {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
    char *const buffer = new char[status._ucount + 1];
    MPI_Recv(buffer, status._ucount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &status);
    buffer[status._ucount] = '\0';
    boost::json::object json = boost::json::parse(buffer).as_object();
    delete[] buffer;
    if (json["type"] == "register") {
      json.clear();
      json["type"] = "ack";
      ++num_registered;
      std::stringstream s;
      s << json;
      auto raw = s.str();
      MPI_Send(raw.c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
    } else if (json["type"] == "deregister") {
      --num_registered;
    } else if (json["type"] == "gridRequest") {
      auto json_bin_widths = json.at("spacing").as_array();
      auto json_node_counts = json.at("nodeCounts").as_array();
      boost::json::object json_to_send;
      double start[3];
      double binwidths[3];
      int node_counts[3];
      start[0] = json.at("offset").as_array().at(0).as_double();
      start[1] = json.at("offset").as_array().at(1).as_double();
      start[2] = json.at("offset").as_array().at(2).as_double();
      binwidths[0] = json_bin_widths.at(0).as_double();
      binwidths[1] = json_bin_widths.at(1).as_double();
      binwidths[2] = json_bin_widths.at(2).as_double();
      node_counts[0] = json_node_counts.at(0).as_int64();
      node_counts[1] = json_node_counts.at(1).as_int64();
      node_counts[2] = json_node_counts.at(2).as_int64();
      json_to_send["nodes"] = boost::json::array();

      bool first_flag = true;
      for (uint x_bin = 0; x_bin < node_counts[0]; ++x_bin) {
        for (uint y_bin = 0; y_bin < node_counts[1]; ++y_bin) {
          for (uint z_bin = 0; z_bin < node_counts[2]; ++z_bin) {
            boost::json::object to_append;
            double pos[3];
            pos[0] = start[0] + binwidths[0] * x_bin;
            pos[1] = start[1] + binwidths[1] * y_bin;
            pos[2] = start[2] + binwidths[2] * z_bin;
            to_append["xIndex"] = x_bin;
            to_append["yIndex"] = y_bin;
            to_append["zIndex"] = z_bin;
            double forces[3];
            forces[0] = forces[1] = forces[2] = 0.0;

            _get_forces(json["atoms"], first_flag, pos, forces);
            first_flag = false;

            to_append["dfx"] = forces[0];
            to_append["dfy"] = forces[1];
            to_append["dfz"] = forces[2];
            json_to_send.at("nodes").as_array().push_back(to_append);
          }
        }
      }
      std::stringstream s;
      s << json_to_send;
      auto raw = s.str();
      MPI_Send(raw.c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
    }
  } while (num_registered != 0);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&comm);
  MPI_Comm_free(&junk_comm);
  MPI_Finalize();
}
