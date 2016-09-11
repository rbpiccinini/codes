#!/bin/sh
rm -rf VTK*
scp -i ~/master/amazon/piccinini.pem ubuntu@ec2-79-125-86-170.eu-west-1.compute.amazonaws.com:/home/ubuntu/master/sydney/case/VTK.tar.gz .
tar xvfz VTK.tar.gz
