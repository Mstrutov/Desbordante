import desbordante
import pandas
from tabulate import tabulate

TABLE = 'examples/datasets/shipping.csv'

def print_data_frame(data_frame, title = None):
    print_table(data_frame, 'keys', title)

def print_table(table, headers = None, title = None):
    if title is not None:
        print(title)

    print(tabulate(table, headers=headers, tablefmt='psql'))

def print_named_ods(list_ods, data_frame):
    for list_od in list_ods:
        lhs_names = []
        rhs_names = []
        for col_i in list_od.lhs:
            lhs_names.append(data_frame.columns[col_i])
        for col_i in list_od.rhs:
            rhs_names.append(data_frame.columns[col_i])
        print(lhs_names, '->', rhs_names)



if __name__ == '__main__':
    algo = desbordante.od.algorithms.Order()
    algo.load_data(table=(TABLE, ',', True))
    algo.execute()

    ods = algo.get_list_ods()

    df = pandas.read_csv(TABLE)

    print()
    print_data_frame(df)
    print()
    print_named_ods(ods, df)
    print()
    print("Depenency [weight] -> [shipping cost] means that ordering table by weight")
    print("will also order table by shipping cost automatically. Let's order by weight: ")
    
    df_sorted = df.sort_values("weight")
    print()
    print_data_frame(df_sorted)
    print()
    print("We can see that shipping cost is sorted too. And dependency seems reasonable:")
    print("the more the package weighs, the more expensive it will be to send it.")
