.PHONY: ws_code

ws_code:
	@echo "Transfer to .code/Makefile"
	$(MAKE) -f .code/Makefile

clean:
	rm -f *.log *.aux *.md *.out texput.log
