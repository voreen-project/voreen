#copy VE related files

# data dir
INSTALL(FILES  
            data/voreenve.cfg.sample
        DESTINATION data
)

# further directories
INSTALL(DIRECTORY
            resource/voreenve
        DESTINATION resource
)
