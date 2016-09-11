#!/bin/sh
rm -rf VTK*
HOST=ec2-46-51-159-229.eu-west-1.compute.amazonaws.com
scp -i ~/master/amazon/piccinini.pem ubuntu@$HOST:/home/ubuntu/master/sydney/case.coarse/VTK.tar.gz .
tar xvfz VTK.tar.gz
