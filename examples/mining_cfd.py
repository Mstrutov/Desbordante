import desbordante
import pandas

# CFD Discovery parameters:
TABLEPATH = 'examples/datasets/university_fd.csv'
MINIMUM_SUPPORT = 8
MINIMUM_CONFIDENCE = 0.7
MAXIMUM_LHS_COUNT = 3


GREEN_BG_CODE = '\033[1;42m'
RED_BG_CODE = '\033[1;41m'
DEFAULT_BG_CODE = '\033[1;49m'


def items_to_indices(items):
    return [item.attribute for item in items]


def find_var_attributes(items):
    result = []
    for item in items:
        if item.value is None:
            result.append(item.attribute)
    return result


def row_is_supported(row_dataframe, pattern):
    for item in pattern:
        if item.value is not None and str(row_dataframe.iloc[0, item.attribute]).lower() != item.value.lower():
            return False
    return True


# returns a dictionary
#   key: a string representation of a dataframe with certain attribute values
#   value: indices of all rows with those attribute values
def make_lhs_to_row_nums(cfd_lhs_var_attributes, supported_rows, table):
    lhs_to_row_nums = {}
    for row_num in range(0, len(table.index)):
        row_dataframe = table.iloc[[row_num], :]
        if row_num not in supported_rows:
            continue
        if row_dataframe.iloc[0, cfd_lhs_var_attributes].to_string() in lhs_to_row_nums:
            lhs_to_row_nums[row_dataframe.iloc[0, cfd_lhs_var_attributes].to_string()].append(row_num)
        else:
            lhs_to_row_nums[row_dataframe.iloc[0, cfd_lhs_var_attributes].to_string()] = [row_num]
    return lhs_to_row_nums


def get_supported(lhs, table):
    supported = []
    for row_num in range(0, len(table.index)):
        row_dataframe = table.iloc[[row_num], :]
        if row_is_supported(row_dataframe, lhs):
            supported.append(row_num)
            continue
    return supported


def validate_cfd(lhs, rhs, table):
    cfd_lhs_var_attributes = find_var_attributes(lhs)
    supported_rows = get_supported(lhs, table)
    rows_satisfying_cfd = []
    lhs_to_row_nums = make_lhs_to_row_nums(cfd_lhs_var_attributes, supported_rows, table)

    # mark rows with identical lhs with green color if they have the most frequent rhs
    for string_lhs in lhs_to_row_nums:
        # fill rhs_to_row_nums
        rhs_to_row_nums = {}
        for row_num in lhs_to_row_nums[string_lhs]:
            row_dataframe = table.iloc[[row_num], :]
            if row_dataframe.iloc[0, [rhs.attribute]].to_string() in rhs_to_row_nums:
                rhs_to_row_nums[row_dataframe.iloc[0, [rhs.attribute]].to_string()].append(row_num)
            else:
                rhs_to_row_nums[row_dataframe.iloc[0, [rhs.attribute]].to_string()] = [row_num]

        best_count = 0
        most_frequent_rhs = None
        for string_rhs in rhs_to_row_nums:
            if len(rhs_to_row_nums[string_rhs]) > best_count:
                best_count = len(rhs_to_row_nums[string_rhs])
                most_frequent_rhs = string_rhs
        for row_num in rhs_to_row_nums[most_frequent_rhs]:
            rows_satisfying_cfd.append(row_num)

    return rows_satisfying_cfd, supported_rows


def visualize_cfd(cfd, table):
    print('CFD:')
    print(cfd, ":\n")

    rows_satisfying_cfd, supported_rows = validate_cfd(cfd.lhs_items, cfd.rhs_item, table)

    header, *rows = table.to_string().splitlines()
    print(DEFAULT_BG_CODE, header)
    for index, row in enumerate(rows):
        if index not in supported_rows:
            color_code = DEFAULT_BG_CODE
        elif index in rows_satisfying_cfd:
            color_code = GREEN_BG_CODE
        else:
            color_code = RED_BG_CODE
        print(color_code, row, DEFAULT_BG_CODE)

    support = len(supported_rows)
    num_rows_satisfy_cfd = len(rows_satisfying_cfd)
    print("lhs count: ", len(cfd.lhs_items))
    print("support: ", support, "\033[1;42m \033[1;41m \033[1;49m")
    print("confidence: \033[1;32m", num_rows_satisfy_cfd, "\033[1;37m/", support,
          " = ", format(num_rows_satisfy_cfd/support, '.4f'))


if __name__ == '__main__':
    tableDF = pandas.read_csv(TABLEPATH)
    algo = desbordante.cfd.algorithms.Default()
    algo.load_data(table=tableDF)
    algo.execute(cfd_minconf=MINIMUM_CONFIDENCE, cfd_minsup=MINIMUM_SUPPORT, cfd_max_lhs=MAXIMUM_LHS_COUNT)
    result = algo.get_cfds()
    print("options: \nMINIMUM SUPPORT =", MINIMUM_SUPPORT,
          ", MINIMUM CONFIDENCE =", MINIMUM_CONFIDENCE,
          ", MAXIMUM LHS COUNT =", MAXIMUM_LHS_COUNT)
    print("displaying the first five discovered CFD's (or less):\n")
    for cfd in result[:5]:
        visualize_cfd(cfd, tableDF)
        print("\n\n")
