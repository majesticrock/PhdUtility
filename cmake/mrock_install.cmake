include_guard()

macro(mrock_install comp)
    set(_mrock_install_root "${CMAKE_CURRENT_SOURCE_DIR}")
    set(LIBRARY_DEST ${CMAKE_INSTALL_LIBDIR})
    set(INCLUDE_DEST ${CMAKE_INSTALL_INCLUDEDIR})
    set(CONFIG_DEST ${CMAKE_INSTALL_LIBDIR}/cmake/${comp})

    install(TARGETS ${comp}
        EXPORT ${comp}-targets
        ARCHIVE DESTINATION ${LIBRARY_DEST}
        LIBRARY DESTINATION ${LIBRARY_DEST}
    )
    install(DIRECTORY "${_mrock_install_root}/include/" DESTINATION ${INCLUDE_DEST})

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