# from selenium import webdriver
# import time
# driver = webdriver.Edge()
# driver.get('https://www.baidu.com')
# driver.find_element('id','kw').click()
# driver.find_element('id','kw').clear()
# driver.find_element('id','kw').send_keys('webdriver')
# driver.find_element('id','su').click()
# time.sleep(2)
# driver.close()

from selenium import webdriver
import time
driver = webdriver.Edge()
driver.get('https://www.baidu.com')
driver.set_window_size(300, 400)
driver.maximize_window()
print(driver.page_source)
print(driver.title)
time.sleep(3)
driver.refresh()
print(driver.current_url)
driver.get_screenshot_as_file('baidu.png')
# driver.close()
# driver.quit()
