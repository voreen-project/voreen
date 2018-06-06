# copy qt related files

# resources
INSTALL(DIRECTORY
            resource/voreenqt
        DESTINATION resource
)

# docs
INSTALL(DIRECTORY
            doc
        DESTINATION .
        PATTERN "doc/docbook" EXCLUDE
)