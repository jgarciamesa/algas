# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmai.com>

libsasi_sources = files([
	'fasta.cpp',
	'gap.cpp',
	'utils.cpp',
	'sequence.cpp'
])

libsasi_deps = [cli_dep, doctest_dep]

libsasi = static_library('libsasi', [libsasi_sources],
	include_directories : inc,
	dependencies : libsasi_deps,
	cpp_args : ['-DDOCTEST_CONFIG_DISABLE']
)
