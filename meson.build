project('nanotpt', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14', 'buildtype=release'])

pthread_dep = dependency('threads')

executable('nanotpt',
           sources: ['src/nanotpt.cpp'],
           dependencies: pthread_dep,
           install : false)
