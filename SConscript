Import('env')

import SCons.Node.FS

def create_objs(srcs):
    return [env.Object(src) for src in srcs]

def GlobRecursive(pattern, node='.'):
    results = []
    for f in Glob(str(node) + '/*', source=True):
        if type(f) is SCons.Node.FS.Dir:
            results += GlobRecursive(pattern, f)
    results += Glob(str(node) + '/' + pattern, source=True)
    return results

sources = Glob('src/parser/*.cpp') + Glob('src/model/*.cpp') + Glob('src/util/*.cpp') + Glob('src/algorithms/*.cpp')
objects = create_objs(sources)

metanome = env.Program(target = 'metanome', source = objects + ['src/main.cpp'], LIBS = ['stdc++fs'])
Command("config.json", "../config.json", Copy("$TARGET", "$SOURCE"))

#TESTS
libs = ['stdc++fs', 'gtest', 'pthread']

progr = env.Program('tests/tester', Glob('src/tests/*.cpp') + objects, LIBS = libs)

csv_tests = Glob("src/tests/inputData/*.csv")
for csv in csv_tests:
    Command('tests/inputData/' + csv.name, csv.rstr(), Copy("$TARGET", "$SOURCE"))

test_alias = Alias('test', [progr], progr[0].abspath)
AlwaysBuild(progr)
