################################################################################
# LEVEL_SET_PLUGIN
#
# The level set plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(LEVEL_SET_PLUGIN)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	level_set_plugin.cpp)




if(STATIC)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/level_set ${SOURCES})
else(STATIC)
	add_library(LevelSetPlugin SHARED ${SOURCES})
	target_link_libraries (LevelSetPlugin ug4)
endif(STATIC)
