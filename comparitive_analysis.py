import streamlit as st #type: ignore
from typing import Dict, List, Set, Tuple
import pandas as pd #type: ignore
import numpy as np #type: ignore
from toxicity_scraper import ToxicityDatabaseScraper

class ComparativeVenomAnalyzer:
    """
    Analyzes and compares venom compositions across multiple snake species
    Identifies conserved toxins, unique toxins, and functional differences
    """
    
    def __init__(self):
        self.scraper = ToxicityDatabaseScraper()
        self.venom_profiles = self._build_venom_profiles()
    
    def _build_venom_profiles(self) -> Dict:
        """
        Build comprehensive venom profiles for each snake species
        Maps toxins to their functions and characteristics
        """
        return {
            'Black Mamba': {
                'region': 'Sub-Saharan Africa',
                'venom_type': 'Neurotoxic',
                'toxins': {
                    'Three-Finger Toxins': {'percentage': 45, 'function': 'Neuromuscular blockade'},
                    'Serine Proteases': {'percentage': 20, 'function': 'Anticoagulation'},
                    'Phospholipase A2': {'percentage': 15, 'function': 'Anticoagulation + hemolysis'},
                    'Metalloproteinases': {'percentage': 10, 'function': 'Hemorrhage'},
                    'Other': {'percentage': 10, 'function': 'Miscellaneous'}
                },
                'speed': 'Very Fast',
                'victim_behavior': 'Immediate paralysis, respiratory failure'
            },
            'King Cobra': {
                'region': 'Southeast Asia',
                'venom_type': 'Neurotoxic + Hemorrhagic',
                'toxins': {
                    'Three-Finger Toxins': {'percentage': 35, 'function': 'Neuromuscular blockade'},
                    'Serine Proteases': {'percentage': 25, 'function': 'Anticoagulation'},
                    'Metalloproteinases': {'percentage': 20, 'function': 'Hemorrhage + tissue damage'},
                    'Phospholipase A2': {'percentage': 12, 'function': 'Cardiotoxicity'},
                    'Other': {'percentage': 8, 'function': 'Miscellaneous'}
                },
                'speed': 'Fast',
                'victim_behavior': 'Paralysis, cardiovascular collapse'
            },
            'Inland Taipan': {
                'region': 'Australia',
                'venom_type': 'Neurotoxic + Coagulotoxic',
                'toxins': {
                    'Three-Finger Toxins': {'percentage': 40, 'function': 'Neuromuscular blockade'},
                    'Phospholipase A2': {'percentage': 30, 'function': 'Anticoagulation + paralysis'},
                    'Serine Proteases': {'percentage': 15, 'function': 'Blood clotting disruption'},
                    'Metalloproteinases': {'percentage': 10, 'function': 'Mild hemorrhage'},
                    'Other': {'percentage': 5, 'function': 'Miscellaneous'}
                },
                'speed': 'Very Fast',
                'victim_behavior': 'Paralysis, coagulopathy, organ failure'
            },
            'Indian Cobra': {
                'region': 'South Asia',
                'venom_type': 'Neurotoxic',
                'toxins': {
                    'Three-Finger Toxins': {'percentage': 50, 'function': 'Neuromuscular blockade'},
                    'Serine Proteases': {'percentage': 20, 'function': 'Coagulation disruption'},
                    'Phospholipase A2': {'percentage': 15, 'function': 'Cardiotoxicity'},
                    'Metalloproteinases': {'percentage': 10, 'function': 'Tissue damage'},
                    'Other': {'percentage': 5, 'function': 'Miscellaneous'}
                },
                'speed': 'Moderate',
                'victim_behavior': 'Progressive paralysis'
            },
            'Saw-scaled Viper': {
                'region': 'Africa, Middle East, Asia',
                'venom_type': 'Hemorrhagic + Coagulotoxic',
                'toxins': {
                    'Metalloproteinases': {'percentage': 40, 'function': 'Hemorrhage, hemorrhagic shock'},
                    'Serine Proteases': {'percentage': 30, 'function': 'Coagulopathy, hemolysis'},
                    'Phospholipase A2': {'percentage': 15, 'function': 'Hemolysis, kidney damage'},
                    'Three-Finger Toxins': {'percentage': 10, 'function': 'Mild neurotoxicity'},
                    'Other': {'percentage': 5, 'function': 'Miscellaneous'}
                },
                'speed': 'Moderate',
                'victim_behavior': 'Hemorrhage, organ failure'
            },
            'Puff Adder': {
                'region': 'Sub-Saharan Africa',
                'venom_type': 'Hemorrhagic + Cytotoxic',
                'toxins': {
                    'Metalloproteinases': {'percentage': 45, 'function': 'Hemorrhage, tissue destruction'},
                    'Serine Proteases': {'percentage': 25, 'function': 'Coagulation disruption'},
                    'Phospholipase A2': {'percentage': 15, 'function': 'Hemolysis, myotoxicity'},
                    'Three-Finger Toxins': {'percentage': 5, 'function': 'Weak neurotoxicity'},
                    'Other': {'percentage': 10, 'function': 'Miscellaneous'}
                },
                'speed': 'Moderate',
                'victim_behavior': 'Severe tissue damage, necrosis'
            },
            'Rattlesnake': {
                'region': 'North America',
                'venom_type': 'Hemorrhagic + Neurotoxic',
                'toxins': {
                    'Metalloproteinases': {'percentage': 40, 'function': 'Hemorrhage, tissue damage'},
                    'Serine Proteases': {'percentage': 30, 'function': 'Anticoagulation, hemolysis'},
                    'Phospholipase A2': {'percentage': 15, 'function': 'Hemolysis, myotoxicity'},
                    'Neurotoxins': {'percentage': 10, 'function': 'Mild neurotoxicity'},
                    'Other': {'percentage': 5, 'function': 'Miscellaneous'}
                },
                'speed': 'Moderate',
                'victim_behavior': 'Hemorrhage, some paralysis'
            }
        }
    
    def get_all_snakes(self) -> List[str]:
        """Get list of all available snakes"""
        return list(self.venom_profiles.keys())
    
    def compare_toxins(self, snakes: List[str]) -> Dict:
        """
        Compare toxins across multiple snakes
        Returns: unique, conserved, and shared toxins
        """
        if not snakes or len(snakes) < 2:
            return {'error': 'Select at least 2 snakes to compare'}
        
        # Get all toxin types for each snake
        snake_toxins = {}
        for snake in snakes:
            if snake in self.venom_profiles:
                toxins = set(self.venom_profiles[snake]['toxins'].keys())
                snake_toxins[snake] = toxins
        
        # Find conserved toxins (present in all selected snakes)
        if snake_toxins:
            conserved = set.intersection(*snake_toxins.values())
        else:
            conserved = set()
        
        # Find unique toxins for each snake
        all_toxins = set()
        for toxins in snake_toxins.values():
            all_toxins.update(toxins)
        
        unique_by_snake = {}
        for snake in snakes:
            if snake in snake_toxins:
                unique = snake_toxins[snake] - conserved
                unique_by_snake[snake] = unique
        
        return {
            'conserved': conserved,
            'unique_by_snake': unique_by_snake,
            'all_toxins': all_toxins,
            'snake_toxins': snake_toxins
        }
    
    def create_comparison_dataframe(self, snakes: List[str]) -> pd.DataFrame:
        """
        Create a DataFrame comparing key metrics across snakes
        """
        comparison_data = []
        
        for snake in snakes:
            if snake in self.venom_profiles:
                profile = self.venom_profiles[snake]
                
                # Get clinical data
                clinical_db = self.scraper.get_clinical_effects_database()
                clinical = clinical_db.get(snake, {})
                
                comparison_data.append({
                    'Snake Species': snake,
                    'Region': profile['region'],
                    'Venom Type': profile['venom_type'],
                    'Speed': profile['speed'],
                    'LD50 (mg/kg)': clinical.get('ld50_mice', 'N/A'),
                    'Mortality': clinical.get('mortality_rate', 'N/A'),
                    'Onset Time': clinical.get('onset_time', 'N/A'),
                    'Main Toxin': max(profile['toxins'].items(), 
                                     key=lambda x: x[1]['percentage'])[0]
                })
        
        return pd.DataFrame(comparison_data)
    
    def create_toxin_percentage_data(self, snakes: List[str]) -> pd.DataFrame:
        """
        Create DataFrame with toxin percentages for visualization
        """
        data = []
        
        for snake in snakes:
            if snake in self.venom_profiles:
                toxins = self.venom_profiles[snake]['toxins']
                for toxin_name, toxin_data in toxins.items():
                    data.append({
                        'Snake': snake,
                        'Toxin Type': toxin_name,
                        'Percentage': toxin_data['percentage']
                    })
        
        return pd.DataFrame(data)
    
    def create_venom_profile_summary(self, snake: str) -> Dict:
        """
        Get detailed profile for a specific snake
        """
        if snake not in self.venom_profiles:
            return None
        
        profile = self.venom_profiles[snake]
        clinical_db = self.scraper.get_clinical_effects_database()
        clinical = clinical_db.get(snake, {})
        
        return {
            'name': snake,
            'profile': profile,
            'clinical': clinical
        }
    
    def get_evolutionary_insights(self, snakes: List[str]) -> Dict:
        """
        Provide evolutionary and functional insights about the selected snakes
        """
        venom_types = set()
        regions = set()
        toxin_strategies = {}
        
        for snake in snakes:
            if snake in self.venom_profiles:
                profile = self.venom_profiles[snake]
                venom_types.add(profile['venom_type'])
                regions.add(profile['region'])
                
                # Identify primary strategy
                primary_toxin = max(profile['toxins'].items(), 
                                   key=lambda x: x[1]['percentage'])
                strategy = primary_toxin[0]
                toxin_strategies[snake] = strategy
        
        insights = {
            'venom_types': venom_types,
            'regions': regions,
            'primary_strategies': toxin_strategies,
            'diversity': len(venom_types)
        }
        
        # Add evolutionary narrative
        if len(snakes) == 2:
            snake1, snake2 = snakes[0], snakes[1]
            profile1 = self.venom_profiles.get(snake1, {})
            profile2 = self.venom_profiles.get(snake2, {})
            
            toxins1 = set(profile1.get('toxins', {}).keys())
            toxins2 = set(profile2.get('toxins', {}).keys())
            shared = toxins1.intersection(toxins2)
            
            insights['comparison'] = {
                'snake1': snake1,
                'snake2': snake2,
                'shared_toxins': len(shared),
                'total_shared': shared,
                'divergence': len(toxins1 - toxins2) + len(toxins2 - toxins1)
            }
        
        return insights
    
    def get_toxin_function_analysis(self, snakes: List[str]) -> Dict:
        """
        Analyze the functional categories of toxins in selected snakes
        """
        function_analysis = {
            'Neuromuscular Blockade': [],
            'Hemorrhage': [],
            'Coagulopathy': [],
            'Hemolysis': [],
            'Cardiotoxicity': [],
            'Tissue Damage': [],
        }
        
        for snake in snakes:
            if snake in self.venom_profiles:
                profile = self.venom_profiles[snake]
                for toxin_name, toxin_data in profile['toxins'].items():
                    function = toxin_data['function']
                    
                    if 'Neuromuscular' in function or 'paralysis' in function.lower():
                        function_analysis['Neuromuscular Blockade'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
                    if 'Hemorrhage' in function or 'hemorrhagic' in function.lower():
                        function_analysis['Hemorrhage'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
                    if 'Coagulation' in function or 'Coagulopathy' in function:
                        function_analysis['Coagulopathy'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
                    if 'Hemolysis' in function or 'hemolysis' in function.lower():
                        function_analysis['Hemolysis'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
                    if 'Cardiotoxicity' in function or 'Cardiac' in function:
                        function_analysis['Cardiotoxicity'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
                    if 'Tissue' in function or 'necrosis' in function.lower():
                        function_analysis['Tissue Damage'].append(
                            {'snake': snake, 'toxin': toxin_name, 'percentage': toxin_data['percentage']}
                        )
        
        return function_analysis


@st.cache_data(ttl=86400)
def get_comparative_analyzer():
    """
    Cached instance of the analyzer
    """
    return ComparativeVenomAnalyzer()