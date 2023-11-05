# RRTMG_GridComp/CMakeLists.txt

add_library(RRTMG_GridComp STATIC
    RRTMG_GridCompMod.F90
)
add_subdirectory(RRTMG)
target_link_libraries(RRTMG_GridComp PUBLIC GeosCore PRIVATE RRTMGCore)
