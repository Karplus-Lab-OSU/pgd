import unittest
import time
from selenium import webdriver


class SearchPage(unittest.TestCase):
    def setUp(self):

       # Create a new instance of the Firefox driver
        self.driver = webdriver.Firefox()

    def test_removed_options_persist(self):

        # Load search page
        self.driver.get("http://localhost:8000/search")

        # Select the box that indicates number of residues
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Select the first amino acid and click it.
        column2 = self.driver.find_element_by_id("id_aa_choices_list_col_2")
        option = column2.find_elements_by_tag_name('li')[0]
        option.click()

        #Hackish way to do it, but there doesn't seem to be any other
        #common ways to do it.
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "3":
                option.click()

        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        column2 = self.driver.find_element_by_id("id_aa_choices_list_col_2")
        option = column2.find_elements_by_tag_name('li')[0]

        #See if the amino acid selection persisted.
        self.assertFalse(option.get_attribute("class") == "selected")


if __name__ == "__main__":
    unittest.main()
