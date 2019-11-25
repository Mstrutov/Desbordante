AddOption('--release', dest='release', nargs=0)

if GetOption('release') == ():
    env = Environment(
        CXX = 'g++',
        CPPFLAGS = ['-std=c++17', '-O3'])
else:
    env = Environment(
        CXX = 'g++',
        CPPFLAGS = [ '-O0', '-std=c++17',
        #'-fvar-tracking',
        '-ggdb', '-fstack-protector-all','-pedantic', '-Wall', '-Werror',
        '-Wextra', '-Wcast-align', '-Wcast-qual',
        #'-Wctor-dtor-privacy', '-Wswitch',
        '-Wdisabled-optimization', '-Wformat=2', '-Winit-self',
        #'-Wlogical-op',
        '-Wmissing-declarations',
        #'-Wnoexcept', '-Wold-style-cast', '-Woverloaded-virtual',
        '-Wredundant-decls', '-Wsign-promo',
        #'-Wstrict-null-sentinel', #'-Wsign-conversion',
        '-Wstrict-overflow=5', '-Wundef', '-Wno-unused', '-Wno-unused-parameter',
        #'-fsanitize=undefined', '-Wno-gnu-zero-variadic-macro-arguments'           # for clang
        ])

env.Append(CPPPATH=['src'] )

SConscript('SConscript', variant_dir = 'build', exports=['env'], duplicate=0)
