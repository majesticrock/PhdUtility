include_guard()

macro(mrock_install comp)
    set(LIBRARY_DEST lib)
    set(INCLUDE_DEST include)
    set(CONFIG_DEST lib/cmake/${comp})

    install(TARGETS ${comp}
        EXPORT ${comp}-targets
        ARCHIVE DESTINATION ${LIBRARY_DEST}
        LIBRARY DESTINATION ${LIBRARY_DEST}
    )
    install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/ DESTINATION ${INCLUDE_DEST})

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        "${CMAKE_CURRENT_LIST_DIR}/cmake/${comp}-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config.cmake"
        INSTALL_DESTINATION ${CONFIG_DEST}
    )
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config-version.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${comp}-config-version.cmake"
        DESTINATION ${CONFIG_DEST}
    )

    install(EXPORT ${comp}-targets
        FILE ${comp}-targets.cmake
        NAMESPACE mrock::
        DESTINATION ${CONFIG_DEST}
    )
endmacro()