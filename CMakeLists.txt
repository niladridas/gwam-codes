cmake_minimum_required(VERSION 2.6)
project(Mycodes)

set(MY_CODES_DIRS /home/niladriisl/eclipse_workspace/gwam-codes)
#set(SAM_LIBS QuadTS constants)
set(SAM_LIBS constants samlibs)





## libbarrett
find_package(Barrett REQUIRED)
include_directories(${BARRETT_INCLUDE_DIRS} ${MY_CODES_DIRS}/include)
link_directories(${BARRETT_LIBRARY_DIRS} ${MY_CODES_DIRS}/builds)
add_definitions(${BARRETT_DEFINITIONS})
add_definitions(-o2 -g)
## curses
find_package(Curses REQUIRED)
include_directories(${CURSES_INCLUDE_DIR})


#add_library(QuadTS ${MY_CODES_DIRS}/samLibs/QuadTS.cpp)
add_library(constants ${MY_CODES_DIRS}/samLibs/constants.cpp)
add_library(samlibs ${MY_CODES_DIRS}/samLibs/samlibs.cpp)
#add_library(Dynamics ${MY_CODES_DIRS}/samLibs/Dynamics.cpp)

function(add_programs)
	foreach (prog ${ARGV})
		add_executable(${prog} src/${prog}.cpp)
		target_link_libraries(${prog} ${BARRETT_LIBRARIES} ${SAM_LIBS})
	endforeach()
endfunction()

add_programs(
#	test_ex04
#	test_Controller
#	test_dataLog
#	dataLog_jaType
#	collectData
#	collectData2
#	computeCost_rt
#	test_dataLog_running
#	test_Controller1
#	test_Controller2
#	test_Controller3
#	playTenisLog
#	collectData3
#	collectData_onlyB
#	wamMoveNonThread
#	test_myController
#	test_Controller_compController
#	hand_test
#	dataCollector_myController
#	gmm_traj_test
#	test_gmm
	torque_aritra_nilu










)

#target_link_libraries(test_ex04 ${CURSES_LIBRARIES})
