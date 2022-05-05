# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmai.com>

inc = include_directories(['.', 'include'])

boost_dep = dependency('boost',
	version : '>=1.47.0',
	required : true
)

subdir('lib')

executable('salsa', '../salsa.cpp',
	include_directories : inc,
	link_with : [libsalsa]
)
