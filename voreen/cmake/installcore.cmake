# copy all core data
# root dir
INSTALL(FILES   
            ../Changelog.txt
            ../CREDITS.txt
            ../LICENSE.txt
            ../LICENSE-academic.txt
        DESTINATION .
)

# tgt
INSTALL(DIRECTORY
            ext/tgt/glsl
        DESTINATION ext/tgt
)

# screenshots
INSTALL(DIRECTORY   
            data/screenshots
        DESTINATION data
        PATTERN "data/screenshots/*" EXCLUDE
)

# resources
INSTALL(DIRECTORY
            resource/voreencore
        DESTINATION resource
)
