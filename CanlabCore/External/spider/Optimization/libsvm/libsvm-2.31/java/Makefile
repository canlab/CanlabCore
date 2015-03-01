FILES =	svm.class svm_model.class svm_node.class svm_parameter.class \
	svm_problem.class svm_train.class svm_predict.class svm_toy.class

#JAVAC = jikes
#JAVAC_FLAGS = -classpath /opt/IBMJava2-13/jre/lib/rt.jar
JAVAC = javac
JAVAC_FLAGS =

all: $(FILES)
	jar cvf libsvm.jar *.class
	
%.class: %.java
	$(JAVAC) $(JAVAC_FLAGS) $<

svm.java: svm.m4
	m4 svm.m4 > svm.java

clean:
	rm -f *.class *.jar *~ svm.java

dist: clean all
	rm *.class

	