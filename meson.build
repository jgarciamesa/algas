# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>

# 	ALignment GAp Statistic
# Init project
project('sasi', 'cpp',
	version : '0.1.9000',
	license : 'MIT',
	meson_version : '>=0.58.0',
	default_options : ['buildtype=debugoptimized', 'cpp_std=c++17'])

subdir('contrib')
subdir('src')
subdir('tests')

configure_file(input : 'sasi.sublime-project.in',
	output : 'sasi.sublime-project',
	configuration : {
		'source_dir' : meson.project_source_root(),
		'project_name' : meson.project_name(),
		'build_dir' : meson.project_build_root()
	}
)