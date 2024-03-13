import click

class TableParamType(click.ParamType):
    name = "TABLE"

    def convert(self, value, param, ctx) \
            -> tuple[str, str, bool]:
        if isinstance(value, tuple[str, str, bool]):
            return value

        str_tuple = value.split()
        if len(str_tuple) != 3:
            self.fail(f"ERROR: {value!r} is not a valid table description", param, ctx)
        return (str_tuple[0], str_tuple[1], bool(str_tuple[3]))
