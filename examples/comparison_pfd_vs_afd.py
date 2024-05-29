import desbordante
from ordered_set import OrderedSet

TABLE = 'examples/datasets/glitchy_sensor.csv'
ERROR = 0.18
ERROR_MEASURE = 'per_value' # per_tuple or per_value


def stringify(fds):
    return OrderedSet(map(str, fds))


def get_afds():
    algo = desbordante.afd.algorithms.Tane()
    algo.load_data(table=(TABLE, ',', True))
    algo.execute(error=ERROR)
    return algo.get_fds()


def get_pfds():
    algo = desbordante.pfd.algorithms.PFDTane()
    algo.load_data(table=(TABLE, ',', True))
    algo.execute(error=ERROR, error_measure=ERROR_MEASURE)
    return algo.get_fds()


pfds = OrderedSet(get_pfds())
afds = OrderedSet(get_afds())

print("pFDs \ AFDs =", stringify(pfds - afds))
print("AFDs \ pFDs =", stringify(afds - pfds))
print("AFDs ∩ pFDs =", stringify(afds & pfds))

print("1 - PerValue([DeviceId] -> Data) =", 0.1714285714)

verifier_algo = desbordante.afd_verification.algorithms.Default()
verifier_algo.load_data(table=(TABLE, ',', True))
for fd in pfds - afds:
    verifier_algo.execute(lhs_indices=fd.lhs_indices, rhs_indices=[fd.rhs_index])
    fd_error = verifier_algo.get_error()
    print(f"e({fd}) =", fd_error) # AFD error is signifcantly larger than PFD PerValue

print('In case of PerValue error measure, violations on data from the single "glitchy"')
print('sensor device among many do not prevent dependecy from being found')
