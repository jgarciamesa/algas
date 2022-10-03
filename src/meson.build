# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmai.com>

inc = include_directories(['.', 'include'])

subdir('lib')

executable('sasi', 'sasi.cpp',
	include_directories : inc,
	dependencies : [libsasi_deps],
	link_with : [libsasi],
	install : true
)
