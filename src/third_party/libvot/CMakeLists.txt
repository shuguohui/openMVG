PROJECT (libvot CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

ADD_LIBRARY (libvot clustering.h clustering.cpp data_structures.h  data_structures.cpp data_types.h global_params.h global_params.cpp io_utils.h io_utils.cpp vocab_tree.h vocab_tree.cpp vot_pipeline.h vot_pipeline.cpp)
SET_PROPERTY(TARGET libvot PROPERTY FOLDER OpenMVG/3rdParty)

if(LIBVOT_INCLUDE_INSTALL_DIR)
  set(INCLUDE_INSTALL_DIR ${LIBVOT_INCLUDE_INSTALL_DIR})
else()
  set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/libvot")
endif()

INSTALL(
  DIRECTORY .
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h"
)
INSTALL(TARGETS libvot DESTINATION lib EXPORT openMVG-targets)

