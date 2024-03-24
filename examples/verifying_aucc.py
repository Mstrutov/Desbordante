import desbordante
import pandas as pd
import termcolor

def print_results(ucc_verifier):
    if ucc_verifier.ucc_holds():
        print('UCC holds, showing stats for AUCC is useless')
    else:
        print('UCC does not hold')
        print(f'But AUCC with error = {"{:.4f}".format(algo.get_error())} holds')
        print()
        print('Also:')
        print(f'Total number of rows violating UCC: {ucc_verifier.get_num_rows_violating_ucc()}')
        print(f'Number of clusters violating UCC: {ucc_verifier.get_num_clusters_violating_ucc()}')
        print('Clusters violating UCC:')
        print(f'found {ucc_verifier.get_num_clusters_violating_ucc()} clusters violating UCC:')
        table = pd.read_csv('examples/datasets/AUCC_example.csv')
        columns = []
        violating_rows = sum(ucc_verifier.get_clusters_violating_ucc(),start=[])
        for i, row in table.iterrows():
            if i not in violating_rows:
                print(*list(row), sep="\t")
            else:
                print(
                    *list(
                        map(
                            lambda x: termcolor.colored(str(x), "red"),
                            row,
                        )
                    ), sep="\t"
                )
        # НАЙДЕНО X кластеров нарушающих UCC, для каждого кластера достать все строки, подсветить красным значения
        # Печатать только найденные кластера
        clusters_violating_ucc = ucc_verifier.get_clusters_violating_ucc()
        for cluster in clusters_violating_ucc:
            print(cluster)        
        
    print()


# Loading input data
algo = desbordante.aucc_verification.algorithms.Default()
algo.load_data(table=('examples/datasets/AUCC_example.csv', ',', True))

# Checking whether (ID) UCC holds
algo.execute(ucc_indices=[0])
print('-'*80)
print('Checking whether (ID) UCC holds')
print('-'*80)
print()
print_results(algo)

# Checking whether (name) UCC holds
algo.execute(ucc_indices=[1])
print('-'*80)
print('Checking whether (name) UCC holds')
print('It should not hold, there are 2 person, named Alex')
print('-'*80)
print()
print_results(algo)

# Checking whether (card_num) UCC holds
algo.execute(ucc_indices=[2])
print('-'*80)
print('Checking whether (card_num) UCC holds')
print('It should not hold, there are 2 same card numbers')
print('-'*80)
print()
print_results(algo)

# Checking whether (card_num, card_active) UCC holds
algo.execute(ucc_indices=[2,3])
print('-'*80)
print('Checking whether (card_num, card_active) UCC holds')
print('It should hold, same cards are not active simultaneously')
print('-'*80)
print()
print_results(algo)
print('-'*80)
