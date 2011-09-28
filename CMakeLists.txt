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



add_library(LevelSetPlugin SHARED ${SOURCES})

if(STATIC)
	target_link_libraries (LevelSetPlugin ug4_s)
else(STATIC)
	target_link_libraries (LevelSetPlugin ug4)
endif(STATIC)
