#+--------------------------------------------------+
#|  R matching with spluss NULL values              |
#|                                                  |
setClassUnion("functionorNULL", c("function", "NULL"))
setClassUnion("expressionorNULL", c("expression", "NULL"))
setClassUnion("callorNULL", c("formula","call", "NULL"))
setClassUnion("integerorNULL", c("integer", "NULL"))
setClassUnion("numericorNULL", c("numeric", "NULL"))
setClassUnion("characterorNULL", c("character", "NULL"))
setClassUnion("logicalorNULL", c("logical", "NULL"))
setClassUnion("listorNULL", c("list", "NULL"))
setClassUnion("nl.formorNULL", c("nl.form", "NULL"))
setClassUnion("nl.fittorNULL", c("nl.fitt", "NULL"))
setClassUnion("nl.fitt.roborNULL", c("nl.fitt.rob", "NULL"))
setClassUnion("fittmethodorNULL", c("fittmethod", "NULL"))