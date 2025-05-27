/**
 * @brief Defines library functions for use in the ARBFN library
 * @author J Dehmel, J Schiffbauer, 2024, MIT License
 */

#include "interchange.h"
#include <boost/json/array.hpp>
#include <boost/json/src.hpp>
#include <cstdint>
#include <iostream>
#include <mpi.h>
#include <thread>

/**
 * @brief Turn a JSON object into a std::string
 * @param _what The JSON to stringify
 * @return The string version
 */
std::string json_to_str(boost::json::value _what)
{
  std::stringstream s;
  s << _what;
  return s.str();
}

/**
 * @brief Yields a string JSON version of the given atom
 * @param _what The atom to JSON-ify
 * @return The serialized version of the atom
 */
boost::json::object to_json(const AtomData &_what)
{
  boost::json::object j;

  j["x"] = _what.x;
  j["vx"] = _what.vx;
  j["fx"] = _what.fx;
  j["y"] = _what.y;
  j["vy"] = _what.vy;
  j["fy"] = _what.fy;
  j["z"] = _what.z;
  j["vz"] = _what.vz;
  j["fz"] = _what.fz;

  if (_what.is_dipole) {
    j["mux"] = _what.mux;
    j["muy"] = _what.muy;
    j["muz"] = _what.muz;
  }

  return j;
}

/**
 * @brief Parses some JSON object into raw fix data.
 * @param _to_parse The JSON object to load from
 * @return The deserialized version of the object
 */
FixData from_json(const boost::json::value &_to_parse)
{
  FixData f;

  f.dfx = _to_parse.at("dfx").as_double();
  f.dfy = _to_parse.at("dfy").as_double();
  f.dfz = _to_parse.at("dfz").as_double();

  return f;
}

/**
 * @brief Await an MPI packet for some amount of time, throwing an error if none arrives.
 * @param _max_ms The max number of milliseconds to wait before error
 * @param _into The `boost::json` to save the packet into
 * @param _received_from Where to save the MPI source of the sender
 * @param _comm The MPI communicator to use
 * @return True on success, false on failure
 */
bool await_packet(const double &_max_ms, boost::json::object &_into, unsigned int &_received_from,
                  MPI_Comm &_comm)
{
  bool got_any_packet;
  std::chrono::high_resolution_clock::time_point send_time, now;
  std::string response;
  uint64_t elapsed_us;
  MPI_Status status;
  int flag;
  char *buffer;

  got_any_packet = false;
  send_time = std::chrono::high_resolution_clock::now();
  while (!got_any_packet) {
    // Check for message recv resolution
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _comm, &flag, &status);
    if (flag && status._ucount > 0) {
      buffer = new char[status._ucount + 1];
      MPI_Recv(buffer, status._ucount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, _comm, &status);
      buffer[status._ucount] = '\0';
      response = buffer;
      delete[] buffer;

      got_any_packet = true;
      _received_from = status.MPI_SOURCE;

      break;
    }

    // Update time elapsed
    now = std::chrono::high_resolution_clock::now();
    elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(now - send_time).count();

    // If it has been too long, indicate error
    if (elapsed_us / 1000.0 > _max_ms && _max_ms > 0.0) {
      // Indicate error
      std::cerr << "Timeout!\n";
      return false;
    }

    // Else, sleep for a bit
    std::this_thread::sleep_for(std::chrono::microseconds(250));
  }

  // Unwrap packet
  _into = boost::json::parse(response).as_object();
  return true;
}

/**
 * @brief Send the given atom data, then receive the given fix data. This is blocking, but does not allow worker-side gridlocks.
 * @param _n The number of atoms/fixes in the arrays.
 * @param _from An array of atom data to send
 * @param _into An array of fix data that was received
 * @param _max_ms The max number of milliseconds to await each response
 * @param _controller_rank The rank of the controller within the provided communicator
 * @param _comm The MPI communicator to use
 * @returns true on success, false on failure
 */
bool interchange(const size_t &_n, const AtomData _from[], FixData _into[], const double &_max_ms,
                 const unsigned int &_controller_rank, MPI_Comm &_comm)
{
  bool got_fix, result;
  boost::json::object json_send, json_recv;
  unsigned int received_from;
  boost::json::array list;
  std::string to_send;

  // Prepare and send the packet
  json_send["type"] = "request";
  json_send["expectResponse"] = _max_ms;
  for (size_t i = 0; i < _n; ++i) { list.push_back(to_json(_from[i])); }
  json_send["atoms"] = list;

  to_send = json_to_str(json_send);
  MPI_Send(to_send.c_str(), to_send.size(), MPI_CHAR, _controller_rank, 0, _comm);

  // Await response
  got_fix = false;
  while (!got_fix) {
    // Await any sort of packet
    result = await_packet(_max_ms, json_recv, received_from, _comm);
    if (!result) {
      std::cerr << "await_packet failed\n";
      return false;
    } else if (received_from != _controller_rank) {
      continue;
    }

    // If "waiting" packet, continue. Else, break.
    if (json_recv.at("type") == "waiting") {
      continue;
    } else {
      if (json_recv["type"] != "response") {
        std::cerr << "Controller sent bad packet w/ type '" << json_recv["type"] << "'\n";
        return false;
      }
      got_fix = true;
      break;
    }
  }

  // Transcribe fix data
  if (json_recv.at("atoms").as_array().size() != _n) {
    std::cerr << "Received malformed fix data from controller: Expected " << _n
              << " atoms, but got " << json_recv.at("atoms").as_array().size() << "\n";
    return false;
  }
  for (size_t i = 0; i < _n; ++i) { _into[i] = from_json(json_recv.at("atoms").as_array().at(i)); }

  return true;
}

/**
 * @brief Sends a registration packet to the controller.
 * @return True on success, false on error.
 */
bool send_registration(unsigned int &_controller_rank, MPI_Comm &_comm)
{
  boost::json::object json;
  std::string to_send;
  int world_size, rank;
  bool result;

  MPI_Comm_rank(_comm, &rank);
  MPI_Comm_size(_comm, &world_size);

  json["type"] = "register";
  to_send = json_to_str(json);

  for (int i = 0; i < world_size; ++i) {
    if (i != rank) { MPI_Send(to_send.c_str(), to_send.size(), MPI_CHAR, i, 0, _comm); }
  }

  json.clear();
  do {
    result = await_packet(10000.0, json, _controller_rank, _comm);
    if (!result) { return false; }
  } while (!json.contains("type") || json.at("type") != "ack");

  return true;
}

/**
 * @brief Sends a deregistration packet to the controller.
 */
void send_deregistration(const int &_controller_rank, MPI_Comm &_comm)
{
  std::string to_send = "{\"type\": \"deregister\"}";
  MPI_Send(to_send.c_str(), to_send.size(), MPI_CHAR, _controller_rank, 0, _comm);
}

std::list<FFieldNodeData> ffield_interchange(const double _start[3], const double _bin_widths[3],
                                             const unsigned int _node_counts[3],
                                             const unsigned int &_controller_rank, MPI_Comm &_comm,
                                             uintmax_t &_every,
                                             const unsigned int &_atoms_to_send_size,
                                             const AtomData _atoms_to_send[])
{
  boost::json::object to_send;

  to_send["type"] = "gridRequest";
  to_send["offset"] = boost::json::array({_start[0], _start[1], _start[2]});
  to_send["spacing"] = boost::json::array({_bin_widths[0], _bin_widths[1], _bin_widths[2]});
  to_send["nodeCounts"] = boost::json::array({_node_counts[0], _node_counts[1], _node_counts[2]});

  // Optional section to send atom information
  if (_atoms_to_send_size > 0) {
    boost::json::array list;
    for (size_t i = 0; i < _atoms_to_send_size; ++i) { list.push_back(to_json(_atoms_to_send[i])); }
    to_send["atoms"] = list;
  }

  std::stringstream to_send_strm;
  to_send_strm << to_send;
  const std::string to_send_string = to_send_strm.str();

  MPI_Send(to_send_string.c_str(), to_send_string.size(), MPI_CHAR, _controller_rank, 0, _comm);

  MPI_Status status;
  MPI_Probe(_controller_rank, 0, _comm, &status);

  char *buffer = new char[status._ucount + 1];
  MPI_Recv(buffer, status._ucount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, _comm, &status);
  buffer[status._ucount] = '\0';

  boost::json::object response = boost::json::parse(buffer).as_object();

  // array of points
  std::list<FFieldNodeData> out;

  const auto points = response.at("nodes").as_array();
  for (const auto &point : points) {
    FFieldNodeData to_add;
    to_add.x_index = point.at("xIndex").as_int64();
    to_add.ybin = point.at("yIndex").as_int64();
    to_add.zbin = point.at("zIndex").as_int64();
    to_add.dfx = point.at("dfx").as_double();
    to_add.dfy = point.at("dfy").as_double();
    to_add.dfz = point.at("dfz").as_double();
    out.push_back(to_add);
  }

  delete[] buffer;

  if (response.contains("every")) {
    // Bonus feature: We can change the interval on the fly
    _every = response.at("every").as_uint64();
  }

  return out;
}
