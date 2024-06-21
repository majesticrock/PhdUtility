BUILD_DIR = build
INSTALL_PREFIX = $(HOME)/usr/local

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

# Run cmake to configure the project with custom install prefix
$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX) ..

# Clean up the build directory
clean:
	@rm -rf $(BUILD_DIR)

# Optionally install the header-only library to the custom directory
install: all
	@$(MAKE) -C $(BUILD_DIR) install

.PHONY: all clean install
