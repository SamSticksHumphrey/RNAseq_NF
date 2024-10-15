# Build Docker image for RNAseq pipeline
## S. Humphrey November 2022

To build the docker run:
`docker build -t rnaseq_awsbatch .`

Then to add this to the Amazon ECR:

### Login
Need to sort MFA:
`aws sts get-session-token --serial-number <YourMFADeviceSerialNumberFoundInSecurityCredentials> --token-code <deviceToken>`

Then add the MFA token to your credientials
`aws configure set aws_session_token <token>`

And log into the ECR
`aws ecr get-login-password --region eu-west-2 | docker login --username AWS --password-stdin <AWS number>.dkr.ecr.eu-west-2.amazonaws.com`

### Tag the rnaseq_awsbatch docker to it's name in the ECR
`docker tag  rnaseq_awsbatch:latest 449547545634.dkr.ecr.eu-west-2.amazonaws.com/rnaseq_awsbatch`

### Push the tagged docker to the ECR
`docker push 449547545634.dkr.ecr.eu-west-2.amazonaws.com/rnaseq_awsbatch`
