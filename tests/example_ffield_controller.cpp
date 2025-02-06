/**
 * @brief A simple controller for demonstrating the arbfn/ffield
 * fix. This is not compiled or run by the default test system
 * for the ARBFN repo, but is provided as a basis for ffield
 * controllers.
 */

#include <boost/json/src.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <sstream>
#include <string>

const static double r = 50.0;    // radius about zero for bulk
const static double s = 1.0;     // steepness of sigmoid
const static double m = 0.1;     // magnitude

constexpr double sigmoid(const double _x)
{
  return 1.0 / (1.0 + exp(-_x));
}

int main()
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
    }

    // Everything before this is just like `fix arbfn`

    // The part that's different from a `fix arbfn` controller
    else if (json["type"] == "gridRequest") {
      boost::json::object json_to_send;
      double start[3];
      double binwidths[3];
      int node_counts[3];

      // Get information from the JSON request
      start[1] = json.at("offset").as_array().at(1).as_double();

      auto json_bin_widths = json.at("spacing").as_array();
      auto json_node_counts = json.at("nodeCounts").as_array();

      binwidths[0] = json_bin_widths.at(0).as_double();
      binwidths[1] = json_bin_widths.at(1).as_double();
      binwidths[2] = json_bin_widths.at(2).as_double();

      node_counts[0] = json_node_counts.at(0).as_int64();
      node_counts[1] = json_node_counts.at(1).as_int64();
      node_counts[2] = json_node_counts.at(2).as_int64();

      // Iterate over all needed nodes
      json_to_send["nodes"] = boost::json::array();
      for (uint x_bin = 0; x_bin < node_counts[0]; ++x_bin) {
        for (uint y_bin = 0; y_bin < node_counts[1]; ++y_bin) {
          for (uint z_bin = 0; z_bin < node_counts[2]; ++z_bin) {
            // Create the node to add to the list of all nodes
            boost::json::object to_append;
            to_append["xIndex"] = x_bin;
            to_append["yIndex"] = y_bin;
            to_append["zIndex"] = z_bin;
            to_append["dfx"] = 0.0;
            to_append["dfz"] = 0.0;

            // Calculate the global y position (can also be done
            // for x and z if wanted)
            double y = start[1] + binwidths[1] * y_bin;

            // Calculate the wall effect force field as a sum
            // of sigmoids
            to_append["dfy"] = m * (sigmoid(s * (-y - r)) - sigmoid(s * (y - r)));

            // Append this node
            json_to_send.at("nodes").as_array().push_back(to_append);
          }
        }
      }

      // Send the calculated nodes
      std::stringstream s;
      s << json_to_send;
      auto raw = s.str();
      MPI_Send(raw.c_str(), raw.size(), MPI_CHAR, status.MPI_SOURCE, 0, comm);
    }

    // Everything after this is just like `fix arbfn`
  } while (num_registered != 0);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm_free(&comm);
  MPI_Comm_free(&junk_comm);
  MPI_Finalize();
  return 0;
}
