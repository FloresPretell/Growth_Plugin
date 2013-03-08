################################################################################
# LEVEL_SET
#
# The level set plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(LEVEL_SET)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	level_set_plugin.cpp)




if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/level_set ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(LevelSet SHARED ${SOURCES})
	target_link_libraries (LevelSet ug4)
endif(buildEmbeddedPlugins)
