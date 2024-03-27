import desbordante
import pandas
from tabulate import tabulate

TABLE = 'examples/datasets/salary.csv'
TIME_LIMIT_SECONDS = 3

def print_data_frame(data_frame, title = None):
    print_table(data_frame, 'keys', title)

def print_table(table, headers = None, title = None):
    if title is not None:
        print(title)
    
    print(tabulate(table, headers=headers, tablefmt='psql'))

def print_attribute_symbols(table):
    print('Attribute symbols:')
    
    counter = 1

    for column in table:
        print(f'{column} -- {counter}')
        counter += 1

def print_desc_ods_with_comments(desc_ods):
    print('descending ods:', len(desc_ods))

def print_asc_ods_with_comments(asc_ods, table):
    print('ascending ods:', len(asc_ods))

    for od in asc_ods:
        print(od)

    print()
    print(f'Dependency "{asc_ods[0]}" means that ordering the table')
    print('inside each equivalence class from "year" by attribute "avg_salary"')
    print('automatically entails ordering by attribute "employee_grade".')

    print()
    print('We have 2 equivalence classes in "year": [2020] and [2021].')
    print('Let\'s split the table into two tables based on these classes.')

    table_part1 = table.iloc[:5,:]
    table_part2 = table.iloc[5:,:]

    print()
    print_data_frame(table_part1, 'Part 1: this part of table corresponds to class [2020]')

    print()
    print('Let\'s sort it by attribute "avg_salary".')

    table_part1_sorted = table_part1.sort_values('avg_salary')

    print()
    print_data_frame(table_part1_sorted, 'Sorted part 1:')

    print()
    print('We can see that this sort entails automatic ordering by')
    print('attribute "employee_grade".')

    print()
    print_data_frame(table_part2, 'Part 2: this part of table corresponds to class [2021]')
    
    print()
    print('Let\'s sort it by attribute "avg_salary".')

    table_part2_sorted = table_part2.sort_values('avg_salary')

    print()
    print_data_frame(table_part2_sorted, 'Sorted part 2:')

    print()
    print('We can see that this sort entails automatic ordering by')
    print('attribute "employee_grade" too.')

    print()
    print(f'Dependency "{asc_ods[1]}" is similar to the first and means that')
    print('ordering the table inside each equivalence class from "year" by')
    print('attribute "employee_grade" automatically entails ordering by')
    print('attribute "avg_salary". This can be seen in the tables above.')

    print()
    print('In other words, these dependencies indicate that the ordering of')
    print('average salary entails an automatic ordering of the employee grade')
    print('and vice versa.')

def print_simple_ods_with_comments(simple_ods, table):
    print('simple ods:', len(simple_ods))

    for od in simple_ods:
        print(od)

    print()
    print('These dependencies mean that inside each equivalence class from')
    print('an attribute from their context the constancy of the attribute')
    print('from the right side of the dependency can be traced.')

    employee_grade_classes = [f'[{i}]' for i in table['employee_grade']]
    employee_grade_classes_str = ', '.join(employee_grade_classes)
    
    print()
    print(f'For example, let\'s look at "{simple_ods[0]}". The context of this')
    print('dependency is attribute "employee_grade". We have 8 equivalence classes')
    print(f'in "employee_grade": {employee_grade_classes_str}.')
    print('Since all the elements of attribute "employee_grade" are different,')
    print('each of these classes contains only one element, so constancy within')
    print('each class occurs automatically.')

    print()
    print('To better understand such dependencies, refer to the second example.')

if __name__ == '__main__':
    algo = desbordante.od.algorithms.Fastod()
    algo.load_data(table=(TABLE, ',', True))
    algo.execute(time_limit=TIME_LIMIT_SECONDS)

    asc_ods = algo.get_asc_ods()
    desc_ods = algo.get_desc_ods()
    simple_ods = algo.get_simple_ods()

    table = pandas.read_csv(TABLE)

    print_data_frame(table)
    print()
    print_attribute_symbols(table)
    print()
    print_desc_ods_with_comments(desc_ods)
    print()
    print_asc_ods_with_comments(asc_ods, table)
    print()
    print_simple_ods_with_comments(simple_ods, table)
