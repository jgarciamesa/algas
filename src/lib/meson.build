# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmai.com>

libalgas_sources = files([
	'fasta.cpp',
	'utils.cpp'
])

libalgas_deps = [boost_dep]

libalgas = static_library('libalgas', [libalgas_sources],
	include_directories : inc,
	dependencies : libalgas_deps,
)
