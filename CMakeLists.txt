# RRTMG_GridComp/CMakeLists.txt

add_library(RRTMG_GridComp STATIC
    RRTMG_GridCompMod.F90
)
target_link_libraries(RRTMG_GridComp PUBLIC GeosCore)
