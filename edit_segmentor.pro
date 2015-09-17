include (../../shared.pri)


HEADERS       = edit_segmentor_factory.h \
                edit_segmentor_recover.h
				 
SOURCES       = edit_segmentor_factory.cpp \
                edit_segmentor_recover.cpp

TARGET        = edit_segmentor

RESOURCES     = edit_segmentor.qrc

FORMS         = ui/edit_segmentor_recover.ui

#QT += network
