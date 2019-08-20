add_library(continuation INTERFACE)

target_include_directories(continuation INTERFACE .)

install(FILES Continuation.H ContinuationDecl.H DESTINATION include)
