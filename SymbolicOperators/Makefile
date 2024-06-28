BUILD_DIR = build
CASCADE_LAKE_BUILD_DIR = build_cluster
INSTALL_PREFIX = $(HOME)/usr/local

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX) ..

# Cluster target to build the project with different compiler flags
cluster: $(CASCADE_LAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CASCADE_LAKE_BUILD_DIR)

$(build_cluster)/Makefile: CMakeLists.txt
	@mkdir -p $(CASCADE_LAKE_BUILD_DIR)
	@cd $(CASCADE_LAKE_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCLUSTER_BUILD=ON ..

clean:
	@rm -rf $(BUILD_DIR)

install: all
	@$(MAKE) -C $(BUILD_DIR) install

.PHONY: all clean install
