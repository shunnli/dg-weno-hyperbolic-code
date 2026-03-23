# check_target(targetname): display target properties

function(check_target _target)
    if(NOT TARGET "${_target}")
        message(FATAL_ERROR "${_target} is not a target.")
    endif()

    get_target_property(target_type ${_target} TYPE)
    string(TOLOWER "${target_type}" target_type)
    get_target_property(target_source_dir ${_target} SOURCE_DIR)

    message(STATUS "---------- <Check Target> ----------")
    message(STATUS "name: ${_target}")
    message(STATUS "type: ${target_type}")
    message(STATUS "location: ${target_source_dir}")

    foreach(prop IN ITEMS
        SOURCES
        INCLUDE_DIRECTORIES
        LINK_DIRECTORIES
        LINK_LIBRARIES
        LINK_OPTIONS
        COMPILE_OPTIONS
        COMPILE_DEFINITIONS
    )
        string(TOLOWER "${prop}" prop_lower)

        get_target_property(value ${_target} ${prop})
        if(NOT value STREQUAL "value-NOTFOUND")
            message(STATUS "${prop_lower}:")
            foreach(item IN LISTS value)
                message(STATUS "  * ${item}")
            endforeach()
        endif()

        get_target_property(interface_value ${_target} INTERFACE_${prop})
        if(NOT interface_value STREQUAL "interface_value-NOTFOUND")
            message(STATUS "${prop_lower}: (interface)")
            foreach(item IN LISTS interface_value)
                message(STATUS "  + ${item}")
            endforeach()
        endif()
    endforeach()

    message(STATUS "------------------------------------")
endfunction()
