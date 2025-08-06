all: docker
.PHONY: all dbuild docker

FILE=bfpwr

## build docker image (requires root access for docker)
dbuild: Dockerfile
	docker build \
    -t $(FILE) .

## run RStudio in docker container
docker: dbuild
	echo "open RStudio Server at http://localhost:8787"
	docker run \
    --rm \
	-ti \
	-e DISABLE_AUTH=true \
	-e ROOT=true \
	-e USERID=$(id -u) \
	-e GROUPID=$(id -g) \
	-p 8787:8787 \
	-v $(CURDIR)/paper:/home/rstudio/paper \
	$(FILE)
