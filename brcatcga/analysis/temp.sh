#!/bin/bash

#curl 'https://api.gdc.cancer.gov/legacy/files/bfd9f7ad-f6bd-42dd-9687-dca1a605d35f?pretty=true'

#curl 'https://api.gdc.cancer.gov/legacy/files/ffe640d4-b79e-44d1-ab76-481bd2620561?pretty=true'

curl 'https://api.gdc.cancer.gov/legacy/files?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22file_id%22%2C%22value%22%3A%5B%22bfd9f7ad-f6bd-42dd-9687-dca1a605d35f%22%2C%22ffe640d4-b79e-44d1-ab76-481bd2620561%22%5D%7D%7D&pretty=true'

