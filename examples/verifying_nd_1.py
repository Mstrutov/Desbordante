import desbordante
import csv

def row_to_padded_string(row, widths):
    return ''.join(field.ljust(width) for width, field in zip(widths, row))

def print_table(filename: str, row_numbers = []):
    with open(filename, newline='') as table:
        rows = list(csv.reader(table, delimiter=','))

    needed_rows = [rows[0]]
    if len(row_numbers) > 0:
        for row_num in row_numbers:
            needed_rows.append(rows[row_num + 1])
    else:
        needed_rows = rows

    column_widths = []
    for col_num in range(len(needed_rows[0])):
        max_len = max(len(row[col_num]) for row in needed_rows)
        column_widths.append(max_len + 3)

    header, *data_rows = needed_rows
    print(row_to_padded_string(header, column_widths))
    print('-' * sum(column_widths))
    print('\n'.join(row_to_padded_string(row, column_widths) for row in data_rows))

EXPIRED_PASSPORT_TABLE = ('examples/datasets/nd_verification_datasets/expired_passport.csv', ',', True)
VALID_PASSPORTS_TABLE = ('examples/datasets/nd_verification_datasets/valid_passports.csv', ',', True)

print('Citizens of Arstozka can have no more than two documents: one Arstozka passport and one exit permit.')
print('The following table contains records of some citizens\' documents:')

print_table('examples/datasets/nd_verification_datasets/expired_passport.csv')
print()

algo = desbordante.nd_verification.algorithms.NDVerifier()
algo.load_data(table=EXPIRED_PASSPORT_TABLE)
algo.execute(lhs_indices=[0], rhs_indices=[1], weight=2)

highlights = algo.highlights

print('We need to validate theese data')
print('Let\'s run ND verification algorithm to check that every citizen has no more than two records:')
print(f'\tND holds: {algo.nd_holds}')
print(f'\tWeight is {algo.real_weight}')
print()

print('Let\'s look at highlights (rows that voilate ND):')
highlight_indices = []
for high in highlights:
    highlight_indices += high.get_occurences_indices()

print_table('examples/datasets/nd_verification_datasets/expired_passport.csv', highlight_indices)
print()
print(f'So, {highlights[0].lhs_value} has {highlights[0].occurences_number} documents')
print('One of them is expired and shouldn\'t be in this table. Let\'s remove this line:')
print_table('examples/datasets/nd_verification_datasets/valid_passports.csv')
print()
print('Let\'s run algirhtm again:')

algo = desbordante.nd_verification.algorithms.NDVerifier()
algo.load_data(table=VALID_PASSPORTS_TABLE)
algo.execute(lhs_indices=[0], rhs_indices=[1], weight=2)

print(f'\tND holds: {algo.nd_holds}')