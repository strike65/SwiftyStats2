.PHONY: doc doc-archive test debug release clean

SWIFT_TARGET := SwiftyStats
DOCS_DIR := docs
DOC_ARCHIVE := $(DOCS_DIR)/$(SWIFT_TARGET).doccarchive
SYMBOL_GRAPH_DIR := .build/symbol-graphs
HOSTING_BASE_PATH := swiftystats

doc:
	# Build static site using the DocC SPM plugin
	swift package \
		--allow-writing-to-directory $(DOCS_DIR) \
		generate-documentation \
		--target $(SWIFT_TARGET) \
		--output-path $(DOCS_DIR) \
		--transform-for-static-hosting \
		--hosting-base-path $(HOSTING_BASE_PATH)

doc-archive:
	# Build a .doccarchive using the DocC SPM plugin
	swift package \
		--allow-writing-to-directory $(DOCS_DIR) \
		generate-documentation \
		--target $(SWIFT_TARGET) \
		--output-path $(DOC_ARCHIVE)

clean-docs:
	rm -rf $(DOC_ARCHIVE) $(SYMBOL_GRAPH_DIR) $(DOCS_DIR)

test:
	swift test

debug:
	swift build

release:
	swift build -c release

clean:
	swift package clean
