# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmai.com>

libsalsa_sources = files([
	'fasta.cpp',
	'gap.cpp',
	'utils.cpp',
	'sequence.cpp'
])

libsalsa_deps = [boost_dep]

libsalsa = static_library('libsalsa', [libsalsa_sources],
	include_directories : inc,
	dependencies : libsalsa_deps,
)
