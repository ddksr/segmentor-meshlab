include (../../shared.pri)


HEADERS       = edit_segmentor_factory.h \
                edit_segmentor_recover.h \
                segRecoverDialog.h \
                segMesh.h \
                libsegmentor/*.h
                 
SOURCES       = edit_segmentor_factory.cpp \
                edit_segmentor_recover.cpp \
                segRecoverDialog.cpp \
                segMesh.cpp \
                libsegmentor/*.cpp \
                libsegmentor/*.c

TARGET        = edit_segmentor

RESOURCES     = edit_segmentor.qrc

FORMS         = ui/edit_segmentor_recover.ui

CONFIG+=debug

#QT += network
