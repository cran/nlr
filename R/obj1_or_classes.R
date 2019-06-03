#+--------------------------------------------------+
#|  R matching with spluss NULL values              |
#|                                                  |

setClassUnion("callorNULL", c("formula","call", "NULL"))
setClassUnion("integerorNULL", c("integer", "NULL"))
setClassUnion("numericorNULL", c("numeric", "NULL"))
setClassUnion("characterorNULL", c("character", "NULL"))
setClassUnion("logicalorNULL", c("logical", "NULL"))
setClassUnion("listorNULL", c("list", "NULL"))

setClassUnion("functionorNULL", c("function", "NULL"))
setClassUnion("expressionorNULL", c("expression", "NULL"))
setClassUnion("vectororNULL", c("vector", "NULL"))
setClassUnion("matrixororNULL", c("matrix", "NULL"))
setClassUnion("vectororMatrix", c("vector", "matrix"))