BUILD_DIR = ../build/SymbolicOperators
CLUSTER_BUILD_DIR = ../build/cluster/SymbolicOperators
INSTALL_PREFIX = $(HOME)/usr/local

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

cluster: $(CLUSTER_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CLUSTER_BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX) -DSYM_OP_TESTS=ON ../../SymbolicOperators

$(CLUSTER_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(CLUSTER_BUILD_DIR)
	@cd $(CLUSTER_BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX) -DCLUSTER_BUILD=ON -DBUILD_TESTS=ON ../../SymbolicOperators

test: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	cd $(BUILD_DIR)/tests && ctest --output-on-failure

clean:
	@rm -rf $(BUILD_DIR) $(CLUSTER_BUILD_DIR)

install: all
	@$(MAKE) -C $(BUILD_DIR) install

install_cluster: cluster
	@$(MAKE) -C $(CLUSTER_BUILD_DIR) install

.PHONY: all clean install test
