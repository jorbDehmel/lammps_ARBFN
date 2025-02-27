'''
Resources for building controllers in python
J Dehmel, '25
'''

import time
import json
from mpi4py import MPI


def independent_controller(single_atom_lambda, max_ms: int = 10_000) -> None:
    '''
    A controller wherein every atom's fix is independent
    of every other atom's. This is much more efficient than a
    dependent controller if dependent control is not needed.

    :param single_atom_lambda: The fix to call on every atom,
        with the atom's data being the JSON first arg and the
        resultant force deltas being saved in the second-fourth
        args.
    '''

    MPI_COMM_WORLD: MPI.Comm = MPI.COMM_WORLD

    junk_comm: MPI.Comm = MPI_COMM_WORLD.Split(0, 0)
    comm: MPI.Comm = MPI_COMM_WORLD.Split(56789, 0)

    print('Started controller.')

    # For as long as there are connections left
    request_instance_counter: int = 0
    requests: int = 0
    num_registered: int = 0
    ms_without_response: int = 0

    while True:
        if ms_without_response >= max_ms:
            print('No response for too long!')
            comm.Abort(10)
            break

        time.sleep(0.1)
        ms_without_response += 100
        if ms_without_response % 1_000 == 0:
            print(
                'Seconds without response: '
                f'{ms_without_response / 1_000.0} / '
                f'{max_ms / 1_000.0}')

        status: MPI.Status = MPI.Status()
        res: bool = comm.Iprobe(
            source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

        if res:
            buff = MPI.buffer.allocate(status.count)
            comm.Recv(buff, status.source, status.tag, status)

            j = json.loads(buff)

            buff.release()

            # Bookkeeping
            if j['type'] == 'register':
                num_registered += 1
                raw: bytes = b'{"type": "ack"}'
                comm.Send(raw, status.source, 0)

            elif j['type'] == 'deregister':
                assert num_registered > 0
                num_registered -= 1

            # Data processing
            elif j['type'] == 'request':
                # Periodically update user
                request_instance_counter += 1
                if request_instance_counter % num_registered == 0:
                    request_instance_counter = 0
                    requests += 1

                    if requests % 100 == 0:
                        print(f'Request #{requests}')

                # Determine fix to send back
                l = []
                for item in j['atoms']:
                    # Processing here
                    dfx, dfy, dfz = single_atom_lambda(item)

                    fix = {
                        'dfx': dfx,
                        'dfy': dfy,
                        'dfz': dfz,
                    }
                    l.append(fix)

                # Properly format the response
                json_to_send = {
                    'type': 'response',
                    'atoms': l,
                }

                s = json.dumps(json_to_send)

                comm.Send(s, raw.size(), status.source, 0)

    print('Halting controller')

    MPI_COMM_WORLD.Barrier()

    comm.Free()
    junk_comm.Free()

    MPI.Finalize()


def ffield_controller(get_forces, max_ms: int = 10_000) -> None:
    '''
    ffield controller for more efficient special cases

    :param get_forces: Maps position (first argument) to forces
    (second argument)
    '''

    MPI_COMM_WORLD: MPI.Comm = MPI.COMM_WORLD

    junk_comm: MPI.Comm = MPI_COMM_WORLD.Split(0, 0)
    comm: MPI.Comm = MPI_COMM_WORLD.Split(56789, 0)

    print('Started controller.')

    # For as long as there are connections left
    num_registered: int = 0
    ms_without_response: int = 0
    while True:
        if ms_without_response >= max_ms:
            print('No response for too long!')
            comm.Abort(10)
            break

        time.sleep(0.1)
        ms_without_response += 100
        if ms_without_response % 1_000 == 0:
            print(
                'Seconds without response: '
                f'{ms_without_response / 1_000.0} / '
                f'{max_ms / 1_000.0}')

        status: MPI.Status = MPI.Status()
        res: bool = comm.Iprobe(
            source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

        if res:
            ms_without_response = 0
            buff = MPI.buffer.allocate(status.count)
            comm.Recv(buff, status.source, status.tag, status)

            j = json.loads(buff.tobytes())
            print(f'Controller got probe of type {j["type"]}')

            buff.release()

            # Bookkeeping
            if j['type'] == 'register':
                num_registered += 1
                raw: bytes = b'{"type": "ack"}'
                comm.Send(raw, status.source, 0)

            elif j['type'] == 'deregister':
                assert num_registered > 0
                num_registered -= 1

                if num_registered == 0:
                    break

            elif j['type'] == 'gridRequest':
                json_bin_widths = j['spacing']
                json_node_counts = j['nodeCounts']

                json_to_send = {}
                start = [0] * 3
                binwidths = [0] * 3
                node_counts = [0] * 3

                start[0] = j['offset'][0]
                start[1] = j['offset'][1]
                start[2] = j['offset'][2]
                binwidths[0] = json_bin_widths[0]
                binwidths[1] = json_bin_widths[1]
                binwidths[2] = json_bin_widths[2]
                node_counts[0] = json_node_counts[0]
                node_counts[1] = json_node_counts[1]
                node_counts[2] = json_node_counts[2]

                json_to_send['nodes'] = []

                for x_bin in range(0, node_counts[0]):
                    for y_bin in range(0, node_counts[1]):
                        for z_bin in range(0, node_counts[2]):
                            to_append = {}
                            pos = [0] * 3
                            pos[0] = start[0] + binwidths[0] * x_bin
                            pos[1] = start[1] + binwidths[1] * y_bin
                            pos[2] = start[2] + binwidths[2] * z_bin
                            to_append['xIndex'] = x_bin
                            to_append['yIndex'] = y_bin
                            to_append['zIndex'] = z_bin

                            forces = get_forces(pos)

                            to_append['dfx'] = forces[0]
                            to_append['dfy'] = forces[1]
                            to_append['dfz'] = forces[2]

                            json_to_send['nodes'].append(to_append)

                comm.Send(json.dumps(json_to_send).encode('utf8'),
                          status.source, 0)

    print('Halting controller')

    comm.Free()
    junk_comm.Free()

    MPI.Finalize()
