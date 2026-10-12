
CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -O2

TARGET = ray_trace
TESTA_TARGET = einstein_A
TESTB_TARGET = einstein_B
METRIC_TEST_TARGET = metric_test

SOURCES = ray_trace.cpp
TESTA_SOURCE = tests/einstein_ringA.cpp
TESTB_SOURCE = tests/einstein_ringB.cpp
METRIC_TEST_SOURCE = tests/test_metric.cpp

HEADERS = metric.h rk45_dromandprince.h newton_raphson.h accretion.h rk45_dormandprince2.h

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

$(TESTA_TARGET): $(TESTA_SOURCE) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(TESTA_SOURCE) -o $(TESTA_TARGET)

$(TESTB_TARGET): $(TESTB_SOURCE) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(TESTB_SOURCE) -o $(TESTB_TARGET)

$(METRIC_TEST_TARGET): $(METRIC_TEST_SOURCE) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(METRIC_TEST_SOURCE) -o $(METRIC_TEST_TARGET)

run: $(TARGET)
	./$(TARGET)

test: $(TESTA_TARGET) $(TESTB_TARGET) $(METRIC_TEST_TARGET)
	./$(TESTA_TARGET)
	./$(TESTB_TARGET)
	./$(METRIC_TEST_TARGET)

clean:
	rm -f $(TARGET) $(TESTA_TARGET) $(TESTB_TARGET) $(METRIC_TEST_TARGET)

.PHONY: clean