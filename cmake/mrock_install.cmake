include_guard()

macro(mrock_install comp)
    set(_mrock_install_root "${CMAKE_CURRENT_LIST_DIR}")

    set(CONFIG_DEST
        ${CMAKE_INSTALL_LIBDIR}/cmake/${comp}
    )

    #
    # Install target
    #
    install(
        TARGETS ${comp}
        EXPORT ${comp}-targets
    )

    #
    # Install headers
    #
    if(EXISTS "${_mrock_install_root}/include")
        install(
            DIRECTORY "${_mrock_install_root}/include/"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
    endif()

    #
    # Package config
    #
    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        "${_mrock_install_root}/cmake/${comp}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config.cmake"
        INSTALL_DESTINATION ${CONFIG_DEST}
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config-version.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(
        FILES
            "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config-version.cmake"
        DESTINATION ${CONFIG_DEST}
    )

    #
    # Export targets
    #
    install(
        EXPORT ${comp}-targets
        FILE ${comp}-targets.cmake
        NAMESPACE mrock::
        DESTINATION ${CONFIG_DEST}
    )
endmacro()