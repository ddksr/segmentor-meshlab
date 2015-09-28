include (../../shared.pri)


HEADERS       = edit_segmentor_factory.h \
                edit_segmentor_recover.h \
                segRecoverDialog.h \
                segMesh.h \
                segDrawer.h \
                segMessaging.h \
                segDescriptionsDialog.h \
                libsegmentor/*.h
                 
SOURCES       = edit_segmentor_factory.cpp \
                edit_segmentor_recover.cpp \
                segRecoverDialog.cpp \
                segMesh.cpp \
                segDrawer.cpp \
                segMessaging.cpp \
                segDescriptionsDialog.cpp \
                ../../meshlabplugins/edit_pickpoints/pickedPoints.cpp \
                libsegmentor/*.cpp \
                libsegmentor/*.c

TARGET        = edit_segmentor

RESOURCES     = edit_segmentor.qrc

FORMS         = ui/edit_segmentor_recover.ui \
                ui/edit_segmentor_descriptions.ui


CONFIG+=debug

#QT += network
