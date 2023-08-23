from selenium import webdriver

driver = webdriver.Edge()
driver.get('http://192.168.41.128:8080/mms/login.html')
driver.find_element('xpath','//div[@class="login-top"]/input[1]').send_keys('admin')
driver.find_element('xpath','//div[@class="login-top"]/input[2]').send_keys('admin123')
driver.find_element('xpath','//div[@class="login-top"]/div[@class="forgot"]/input').click()