BUILD_DIR = ../build/utility
INSTALL_PREFIX = $(HOME)/usr/local

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX) -DUTIL_TESTS=ON ../../utility

test: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	cd $(BUILD_DIR)/tests && ctest --output-on-failure

# Clean up the build directory
clean:
	@rm -rf $(BUILD_DIR)

# Optionally install the header-only library to the custom directory
install: all
	@$(MAKE) -C $(BUILD_DIR) install

.PHONY: all clean install test
