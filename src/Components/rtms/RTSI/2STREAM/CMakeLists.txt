esma_set_this()

file (GLOB srcs 2stream_sourcecode/*.[Ffc] 2stream_sourcecode/*.[fF]90)

esma_add_library (${this}
  SRCS ${srcs}
  )

