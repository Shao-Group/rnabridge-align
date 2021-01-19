bin_PROGRAMS = coral

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util
GRAPH_LIB = $(top_builddir)/lib/graph

#coral_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
coral_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
coral_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB)
coral_LDADD = -lgtf -lgraph -lutil

coral_SOURCES = splice_graph.h splice_graph.cc \
				vertex_info.h vertex_info.cc \
				edge_info.h edge_info.cc \
				interval_map.h interval_map.cc \
				config.h config.cc \
				hit.h hit.cc \
				region.h region.cc \
				partial_exon.h partial_exon.cc \
				fragment.h fragment.cc \
				fcluster.h fcluster.cc \
				junction.h junction.cc \
				bundle_base.h bundle_base.cc \
				bundle2.h bundle2.cc \
				path.h path.cc \
				previewer.h previewer.cc \
				bridger.h bridger.cc \
				reference.h reference.cc \
				coral.h coral.cc \
				hyper_node.h hyper_node.cc \
				main.cc
