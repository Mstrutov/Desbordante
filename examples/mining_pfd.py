import desbordante

ERROR_1 = 0.02777777778 # threshold of the Y->X dependency (PerValue metric)
ERROR_2 = 0.09090909091 # threshold of the Y->X dependency (PerTuple metric)
ERROR_3 = 0.273  # threshold of the X->Y dependency (PerTuple metric)  
ERROR_4 = 0.375  # threshold of the X->Y dependency (PerValue metric) 

algo = desbordante.pfd.algorithms.PFDTane()
algo.load_data(table=('examples/datasets/pfd.csv', ',', True))

for ERROR_MEASURE in ['per_value', 'per_tuple']:
    algo.execute(error=ERROR_1, error_measure=ERROR_MEASURE)
    result = algo.get_fds()
    print(ERROR_MEASURE, 'pFDs:')
    for fd in result:
        print(fd)
